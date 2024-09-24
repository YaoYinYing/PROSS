print_banner(){
    echo "-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="
}

print_banner
echo "Checking for advanced Rosetta Flag setting..."
print_banner

flag_adv_flags=''
# attach advanced flag file
if [[ "$adv_flag" != "" ]]; then
    flag_adv_flags="@$(readlink -f $adv_flag) "
fi
print_banner
echo Advanced flags: $flag_adv_flags

print_banner